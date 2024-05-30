
OPTS = /GA /Gy /EHsc /nologo /GS /W4 /std:c++latest /O2 /Oi /Oy /Ob2 /Ox /Ot /MT /D_CRT_SECURE_NO_WARNINGS

OpenBgiLib=BGI\openbgi64.lib

OBJB = BGI\bgi.obj  BGI\client.obj  BGI\graphics.obj \
       BGI\ipc.obj  BGI\server.obj

{bgi\}.c{bgi\}.obj:
    @cl $(OPTS) /c /Fo$@ $<

Polys.exe: Polys.cpp Point.h $(OpenBgiLib)
    @cl $(OPTS) Polys.cpp $(OpenBgiLib)

$(OpenBgiLib):  $(OBJB)
    @lib /NOLOGO /OUT:$@ $(OBJB)
