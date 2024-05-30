#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <bgi/graphics.h>

#include "Point.h"

using namespace std;

constexpr int Gsize = 900;  // Canvas size

std::default_random_engine rgen(std::random_device{}());
std::uniform_int_distribution<int> unidist{1,Gsize - 2};

#pragma warning(disable:4244 4996)

unsigned int scale[] =
{
    0x303398,0x303398,0x303398,0x303398,0x303398,0x303398,0x303398,0x303398,0x303398,
    0x303398,0x4074b7,0x72aed3,0xaad9e9,0xe0f3fa,
    0xfefcbc,0xffe18c,0xffae5a,0xf56d3a,0xdb2d20,
    0xdb2d20,0xdb2d20,0xdb2d20,0xdb2d20,0xdb2d20,0xdb2d20,0xdb2d20
};

int colorScale(double value,
                        size_t count,
                        unsigned int * cols)
{
    double step = 1./(count-1);
    int no = floor(value/step);
    auto cl = [=](int C)
    {
        unsigned int mask[] = { 0xFF0000, 0xFF00, 0xFF };
        unsigned int shft[] = { 16, 8, 0 };
        int c1 = int((cols[no+1]&mask[C])>>shft[C]);
        int c0 = int((cols[no  ]&mask[C])>>shft[C]);
        return (unsigned int)((value - no*step)/step*(c1-c0) + c0)&0xFF;
    };
    return rgb(cl(0),cl(1),cl(2));
}

// To draw polygone
void draw(int color, const std::vector<poly::Point>& a, bool points = true)
{
    if (points)
    {
        setcolor(WHITE);
        for(size_t i = 0; i < a.size(); ++i)
            fillellipse(a[i].x,Gsize-a[i].y,3,3);
    }
    setcolor(color);
    for(size_t i = 0, j = a.size() - 1; i < a.size(); j = i++)
        line(a[j].x,Gsize-a[j].y,a[i].x,Gsize-a[i].y);
}

// To draw kernel
void draw(int color, const std::vector<poly::KernelVertex>& a, bool points = true)
{
    if (points)
    {
        setcolor(WHITE);
        for(size_t i = 0; i < a.size(); ++i)
            fillellipse(a[i].p.x,Gsize-a[i].p.y,3,3);
    }
    setcolor(color);
    for(size_t i = 0, j = a.size() - 1; i < a.size(); j = i++)
        line(a[j].p.x,Gsize-a[j].p.y,a[i].p.x,Gsize-a[i].p.y);
}


int main(int argc, char ** argv)
{
    vector<poly::Point> a;   // Initial set

    int PointCount = 25;
    if (argc > 1) PointCount = atoi(argv[1]);

    if (PointCount == 0) // datafile!!
    {
        // Datafile format:
        //  x_i  y_i           N points coordinates
        ifstream in(argv[1]);
        for(;;)
        {
            double x, y;
            if (in >> x >> y)
            {
                a.emplace_back(poly::Point(x,y));
                PointCount++;
            }
            else break;
        }
    }

    int gd = CUSTOM, gm = CUSTOM_MODE(Gsize,Gsize);;
    initgraph(&gd, &gm, "RGB");

    if (a.size()== 0)  // not filled yet
    for(size_t i = 0; i < PointCount; ++i)
    {
        int x = unidist(rgen), y = unidist(rgen);
        a.emplace_back(poly::Point(x,y));
    }

    if (true) // MonteCarlo
    {
        cleardevice();
        for(auto p: a) fillellipse(p.x,Gsize-p.y,3,3);
        // Build convex hull
        auto hull = poly::convexHull(a);
        // Draw it
        draw(LIGHTMAGENTA,hull);
        setcolor(LIGHTGREEN);

        // Area
        char txt[40];
        sprintf(txt,"S = %.3lf",poly::Area(hull));
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(50,50,txt);

        readkey();

        // MonteCarlo
        std::vector<poly::Point> m, M;
        poly::Point pm, pM;
        double sm = poly::Area(hull), sM = 0;
        for(int i = 0; i < 100000; ++i)
        {
            poly::Point pt(unidist(rgen),unidist(rgen));
            if (!poly::inHull(pt,hull)) { --i; continue; }

            auto r = poly::starPoly(pt,a);
            double s = poly::Area(r);

            if (s < sm)
            {
                cout << (sm = s) << endl;
                m = r;
                pm = pt;
            }
            if (s > sM)
            {
                cout << (sM = s) << endl;
                M = r;
                pM = pt;
            }

        }

        // draw minimal polygone
        setcolor(LIGHTCYAN);
        fillellipse(pm.x,Gsize-pm.y,3,3);
        setcolor(LIGHTCYAN);
        draw(LIGHTCYAN,m);
        sprintf(txt,"S = %.3lf",sm);
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(pm.x,Gsize-pm.y,txt);

        readkey();

        // draw maximal polygone
        cleardevice();
        for(auto p: a) fillellipse(p.x,Gsize-p.y,3,3);
        // Build convex hull
        hull = poly::convexHull(a);
        // Draw it
        draw(LIGHTMAGENTA,hull);
        setcolor(LIGHTGREEN);
        sprintf(txt,"S = %.3lf",poly::Area(hull));
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(50,50,txt);

        setcolor(YELLOW);
        fillellipse(pM.x,Gsize-pM.y,3,3);
        setcolor(YELLOW);
        draw(YELLOW,M);
        sprintf(txt,"S = %.3lf",sM);
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(pM.x,Gsize-pM.y,txt);

        readkey();

    }

    if (true) // kernels exhausting search
    {
        cleardevice();

        auto area = poly::Area(poly::convexHull(a));
        std::vector<poly::Point> m, M;
        double sm = area, sM = 0;

        // First point/kernel
        poly::Point pc;
        for(auto p: a)
        {
            fillellipse(p.x,Gsize-p.y,3,3);
            pc.x += p.x;
            pc.y += p.y;
        }
        pc.x /= a.size();
        pc.y /= a.size();

        auto star = poly::starPoly(pc,a);
        poly::StarKernel krnl(star);


        setfillstyle(SOLID_FILL,BLUE);
        fillellipse(pc.x,Gsize-pc.y,5,5);

        // BFS
        queue<poly::StarKernel> Q;
        unordered_set<poly::StarKernel,poly::StarKernelHash,poly::StarKernelComp> S;
        S.insert(krnl);
        Q.push(krnl);
        while(!Q.empty())
        {
            poly::StarKernel K = Q.front();
            Q.pop();
            // Kernel polygone
            vector<int> d;
            for(auto z: K.k)
            {
                d.push_back(z.p.x);
                d.push_back(Gsize - z.p.y);
            }

            // Save min, max
            double sq = poly::Area(K.s);
            double value = sq/area;
            if (sq < sm)
            {
                cout << (sm = sq) << endl;
                m = K.s;
            }
            if (sq > sM)
            {
                cout << (sM = sq) << endl;
                M = K.s;
            }

            // draw filled kernel
            //setcolor(Color(value));
            setcolor(colorScale(value,size(scale),scale));
            //setfillstyle(SOLID_FILL, Color(value));
            setfillstyle(SOLID_FILL, colorScale(value,size(scale),scale));
            //setfillstyle(SOLID_FILL, rgb(rand() % 255, rand() % 255, rand() % 255));
            fillpoly(int(d.size()/2),d.data());

            for(size_t e = 0; e < K.k.size(); ++e)
            {
                poly::StarKernel K2(K,e);
                if (K2.k.size() == 0) continue;

                if (S.insert(K2).second)
                {
                    //draw(YELLOW,K2.k,false); //Kernel borders?
                    Q.push(K2);
                }
            }
        }

        //for(auto p: a) cout << p.n << " " << p.x << " " << p.y << endl;

        readkey();

        cleardevice();

        char txt[20];

        setcolor(LIGHTCYAN);
        draw(LIGHTCYAN,m);
        sprintf(txt,"S = %.3lf",sm);
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(100,Gsize-100,txt);

        readkey();

        cleardevice();

        setcolor(YELLOW);
        draw(YELLOW,M);
        sprintf(txt,"S = %.3lf",sM);
        settextstyle(SANS_SERIF_FONT,HORIZ_DIR,10);
        outtextxy(100,Gsize-100,txt);

        readkey();

    }

    closegraph();

}

