#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <bit>

namespace poly
{

    double atan3(double y, double x)
    {
        double f = atan2(y,x);
        if (f >= 0) return f;
        return 2*3.141592653589793 + f;
    }


    struct Point
    {
        Point(double x = 0, double y = 0):x(x),y(y),n(N++){}

        double x = 0, y = 0;
        size_t n = 0;
        inline static int N = 0;

        static inline bool less(const Point& a, const Point& b)
        {
            return a.x < b.x || a.x == b.x && a.y < b.y;
        }

        static inline bool cw(const Point& a, const Point& b, const Point& c)
        {
            return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) < 0;
        }

        static inline bool ccw(const Point& a, const Point& b, const Point& c)
        {
            return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) > 0;
        }

        static inline double S0(const Point& a, const Point& b)
        {
            return (a.x*b.y - a.y*b.x)/2.;
        }

    };

    inline bool operator == (const Point& a, const Point& b)
    {
        return a.x == b.x && a.y == b.y;
    }

    // http://e-maxx.ru/algo/convex_hull_graham
    inline std::vector<Point> convexHull(const std::vector<Point>& a)
    {
        std::vector<Point> r;
        if (a.size() <= 1) return r;
        r = a;
        // Sort all points left to right
        sort(r.begin(),r.end(), Point::less);
        Point p1 = r.front(), p2 = r.back();

        std::vector<Point> u{p1} /* Up side */, d{p1} /* Down side */;

        for(size_t i = 1; i < r.size(); ++i)
        {

            if (i == r.size() - 1 || Point::cw(p1, r[i], p2))
            {
                while(u.size()>=2 && !Point::cw(u[u.size()-2],u[u.size()-1],r[i])) u.pop_back();
                u.push_back(r[i]);
            }
            if (i == r.size() - 1 || Point::ccw(p1, r[i], p2))
            {
                while(d.size()>=2 && !Point::ccw(d[d.size()-2],d[d.size()-1],r[i])) d.pop_back();
                d.push_back(r[i]);
            }
        }
        r.clear();
        for(size_t i = 0; i < u.size(); ++i)   r.push_back(u[i]);
        for(size_t i = d.size()-2; i > 0; --i) r.push_back(d[i]);

        return r;
    }

    // All points on the hull
    inline std::vector<Point> convexHullFull(const std::vector<Point>& a)
    {
        std::vector<Point> r;
        if (a.size() <= 1) return r;
        r = a;
        // Sort all points left to right
        sort(r.begin(),r.end(), Point::less);
        Point p1 = r.front(), p2 = r.back();

        std::vector<Point> u{p1} /* Up side */, d{p1} /* Down side */;

        for(size_t i = 1; i < r.size(); ++i)
        {

            if (i == r.size() - 1 || !Point::ccw(p1, r[i], p2))
            {
                while(u.size()>=2 && Point::ccw(u[u.size()-2],u[u.size()-1],r[i])) u.pop_back();
                u.push_back(r[i]);
            }
            if (i == r.size() - 1 || !Point::cw(p1, r[i], p2))
            {
                while(d.size()>=2 && Point::cw(d[d.size()-2],d[d.size()-1],r[i])) d.pop_back();
                d.push_back(r[i]);
            }
        }
        r.clear();
        for(size_t i = 0; i < u.size(); ++i)   r.push_back(u[i]);
        for(size_t i = d.size()-2; i > 0; --i) r.push_back(d[i]);

        return r;
    }

    inline double Area(const std::vector<poly::Point>& a)
    {
        double S = 0;
        for(size_t i = 0; i < a.size(); ++i)
        {
            poly::Point
                p = i ? a[i-1] : a[a.size()-1],
                q = a[i];
            S += (p.x-q.x)*(p.y+q.y);
        }
        return abs(S/2);
    }

    // Build star poligone
    inline std::vector<poly::Point> starPoly(const poly::Point& pt,
                                             const std::vector<poly::Point>& a)
    {
        poly::Point p = pt;
        std::vector<poly::Point> r;
        std::vector<std::tuple<size_t,double,double>> angls;

        // Is this point hull point?
        bool onHull = false;
        auto ch = convexHullFull(a);
        size_t id = 0;
        for(; id < ch.size(); ++id) if (pt == ch[id]) { onHull = true; break; }


        if (onHull)
        {
            // Small shift :)
            id = (id + ch.size()/2)%ch.size();
            p.x += (a[id].x - p.x)*0.001;
            p.y += (a[id].y - p.x)*0.001;
            onHull = false;
        }

        //cout << "Not on hull\n";
        for(size_t i = 0; i < a.size(); ++i)
        {
            angls.emplace_back(i,
                               (a[i].y == p.y && a[i].x == p.x) ? 0 :
                               poly::atan3(a[i].y-p.y,a[i].x-p.x),
                               (a[i].x-p.x)*(a[i].x-p.x)+(a[i].y-p.y)*(a[i].y-p.y));
        }
        sort(angls.begin(),angls.end(),[](const std::tuple<size_t,double,double>& a,
                                          const std::tuple<size_t,double,double>& b)
             {
                 return
                     std::get<1>(a) < std::get<1>(b) ||
                     std::get<1>(a) == std::get<1>(b) && std::get<2>(a) < std::get<2>(b);
             });
        for(size_t i = 0; i < angls.size(); ++i)
            r.emplace_back(a[std::get<0>(angls[i])]);
        return r;
    }

    // Is point in hull
    // Is line from h[0] crossing some segment?
    inline bool inHull(const poly::Point& p, const std::vector<poly::Point>& h)
    {
        if (p.x < h[0].x) return false;
        std::vector<double> angls;
        for(int i = 1; i < h.size(); ++i) angls.push_back(atan2(h[i].y-h[0].y,h[i].x-h[0].x));

        double f = atan2(p.y-h[0].y,p.x-h[0].x);
        if (f > angls[0] || f < angls[angls.size()-1]) return false;

        auto q = std::upper_bound(angls.begin(), angls.end(), f, std::greater<>());

        if (q == angls.begin() || q == angls.end()) return false;

        int b = int(q - angls.begin()) + 1, c = b - 1;

        double s = p.x*(h[b].y - h[c].y) + h[b].x*(h[c].y - p.y) + h[c].x*(p.y - h[b].y);

        return s > 0;

    }

    inline std::vector<poly::Point> BoundRect(const std::vector<poly::Point>& poly)
    {
        auto min_x = poly[0].x;
        auto max_x = poly[0].x;
        auto min_y = poly[0].y;
        auto max_y = poly[0].y;
        for(size_t i = 1; i < poly.size(); ++i)
        {
            min_x = std::min(min_x, poly[i].x);
            max_x = std::max(max_x, poly[i].x);
            min_y = std::min(min_y, poly[i].y);
            max_y = std::max(max_y, poly[i].y);
        }
        return std::vector<poly::Point>{
            {min_x, min_y},
            {min_x, max_y},
            {max_x, max_y},
            {max_x, min_y}
        };
    }

    // Vector product
    inline double leftOf(const poly::Point& p, const poly::Point& begin, const poly::Point& end)
    {
        return (end.x-begin.x)*(p.y-begin.y)-(p.x-begin.x)*(end.y-begin.y);
    }

    inline poly::Point cross(double Lp, const poly::Point& p,
                             double Lq, const poly::Point& q)
    {
        return { (Lp*q.x - Lq*p.x)/(Lp-Lq),(Lp*q.y - Lq*p.y)/(Lp-Lq) };
    }


    struct KernelVertex
    {
        poly::Point p;
        size_t b,e; // begin-end?
    };

    inline std::vector<poly::KernelVertex> buildKernel(const std::vector<poly::Point>& starPoly);


    struct StarKernel // Kernel of star polygone
    {
        std::vector<poly::Point>  s;        // Polygone
        std::vector<poly::KernelVertex> k;  // Kernel for s
        std::vector<std::pair<size_t,size_t>> ixs;
        std::size_t hashvalue;

        StarKernel(const std::vector<poly::Point>& star):s(star),k(buildKernel(star))
        {
            ixs = StarKernelVal();
            hashvalue = StarKernelHash();
        }
        StarKernel(const StarKernel& kernel, size_t edge) // Unfolding of the star edge 
            :s(kernel.s)                                  // recorded at the vertex of the kernel
        {
            //std::cout << "edge = " << edge << "   kernel.k.size() = " << kernel.k.size() << std::endl;
            std::swap(s[kernel.k[edge].b],s[kernel.k[edge].e]);
            k = buildKernel(s);
            ixs = StarKernelVal();
            hashvalue = StarKernelHash();
        }

        std::vector<std::pair<size_t,size_t>> StarKernelVal()
        {
            ixs.clear();
            for(size_t e = 0; e < k.size(); ++e)
            {
                poly::KernelVertex v = k[e];
                ixs.emplace_back(s[v.b].n,s[v.e].n);
            }
            auto it = std::min_element(ixs.begin(),ixs.end());
            std::rotate(ixs.begin(),it,ixs.end());
            return ixs;
        }

        size_t StarKernelHash()
        {
            size_t z = 0;
            for(int i = 0; i < ixs.size(); ++i)
            {
                // There are many various hashes...

                z ^= std::rotr(ixs[i].first,2*i);
                z ^= std::rotr(ixs[i].second,2*i+1);

                // z ^= ixs[i].first  + 0x9e3779b9 + (z<<6) + (z>>2);
                // z ^= ixs[i].second + 0x9e3779b9 + (z<<6) + (z>>2);
            }
            return z;
        }
    };

    struct StarKernelHash
    {
        size_t operator()(const StarKernel& k) const { return k.hashvalue; }
    };

    struct StarKernelComp
    {
        size_t operator()(const StarKernel& k, const StarKernel& m) const
        {
            return k.ixs == m.ixs;
        }
    };


    inline std::vector<poly::KernelVertex> clipBound(const std::vector<poly::KernelVertex>& bound,
                                                     const std::vector<poly::Point>& star,
                                                     size_t b, size_t e)
    {
        std::vector<poly::KernelVertex> result;
        // Go through all the cutoffs counter clockwise
        // p, q - ponts of star poligone
        for(size_t j = bound.size()-1, i = 0; i < bound.size(); j = i++)
        {
            double Lp = leftOf(bound[j].p,star[b],star[e]);
            double Lq = leftOf(bound[i].p,star[b],star[e]);
            // TODO Optimization: Lp = Lq for the next step

            // If the points are on different sides, add the intersection point
            if (Lp < 0 && Lq > 0)
            {
                result.push_back({cross(Lp,bound[j].p,Lq,bound[i].p), b,e});
            }
            if (Lp > 0 && Lq < 0)
            {
                result.push_back({cross(Lp,bound[j].p,Lq,bound[i].p), bound[i].b, bound[i].e});
            }
            if (Lq > 0) result.push_back(bound[i]);
            if (Lq == 0)
                if (Lp > 0)
                    result.push_back(bound[i]);
                else if (Lp < 0)
                    result.push_back({bound[i].p,b,e});
                else
                    if (bound[i].b == (size_t)-1)
                        result.push_back({bound[i].p,b,e});
                    else
                        result.push_back(bound[i]);
        }

        return result;
    }

    inline std::vector<poly::KernelVertex> buildKernel(const std::vector<poly::Point>& starPoly)
    {
        std::vector<poly::Point> rect = BoundRect(starPoly);
        std::vector<poly::KernelVertex> kernel;

        for (const auto& p : rect)
            kernel.push_back({p, (size_t)-1, (size_t)-1});

        for (size_t j = starPoly.size() - 1, i = 0; i < starPoly.size(); j = i++)
        {
            kernel = clipBound(kernel, starPoly, j, i);  // Points numbers!
        }
        return kernel;
    }

}

