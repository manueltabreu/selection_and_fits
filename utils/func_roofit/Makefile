CFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags --glibs` -lRooFit -lRooFitCore -lMinuit 

## https://root-forum.cern.ch/t/compiling-a-roofit-user-pdf-with-g/16266/3

# clang++ -std=c++11 w.cpp -o w `pkg-config --cflags --libs opencv`

EXECUTABLE0 := RooDoubleCBFast
EXECUTABLE1 := GaussExp
EXECUTABLE2 := RooGaussDoubleSidedExp


$(EXECUTABLE0): $(EXECUTABLE0).cc $(EXECUTABLE0).h
	g++ -c RooDoubleCBFast.cc -fPIC $(CFLAGS) $(LDFLAGS)  -I./vdt
	rootcint -f RooDoubleCBFastDict.cxx -c RooDoubleCBFast.h LinkDef.h
	g++ -c RooDoubleCBFastDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)
	gcc -shared -Wl,-soname,libRooDoubleCBFast.so.1 -o libRooDoubleCBFast.so.1.0 RooDoubleCBFast.o RooDoubleCBFastDict.o
	ln -sf libRooDoubleCBFast.so.1.0 libRooDoubleCBFast.so


$(EXECUTABLE1): $(EXECUTABLE1).cxx $(EXECUTABLE1).h
	g++ -c GaussExp.cxx  -fPIC $(CFLAGS) $(LDFLAGS) -lX11
	rootcint -f GaussExpDict.cxx -c GaussExp.h LinkDef.h
	g++ -c GaussExpDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)
	gcc -shared -Wl,-soname,libGaussExp.so.1 -o libGaussExp.so.1.0 GaussExp.o GaussExpDict.o
	ln -sf libGaussExp.so.1.0 libGaussExp.so

$(EXECUTABLE2): $(EXECUTABLE2).cxx $(EXECUTABLE2).h
	g++ -c RooGaussDoubleSidedExp.cxx  -fPIC $(CFLAGS) $(LDFLAGS) -lX11
	rootcint -f RooGaussDoubleSidedExpDict.cxx -c RooGaussDoubleSidedExp.h LinkDef.h
	g++ -c RooGaussDoubleSidedExpDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)
	gcc -shared -Wl,-soname,libRooGaussDoubleSidedExp.so.1 -o libRooGaussDoubleSidedExp.so.1.0 RooGaussDoubleSidedExp.o RooGaussDoubleSidedExpDict.o
	ln -sf libRooGaussDoubleSidedExp.so.1.0 libRooGaussDoubleSidedExp.so






# libRooDoubleCBFast.so : libRooDoubleCBFast.so.1.0
# 	ln -sf libRooDoubleCBFast.so.1.0 libRooDoubleCBFast.so
# 
# libRooDoubleCBFast.so.1.0 : RooDoubleCBFast.o RooDoubleCBFastDict.o
# 	gcc -shared -Wl,-soname,libRooDoubleCBFast.so.1 -o libRooDoubleCBFast.so.1.0 RooDoubleCBFast.o RooDoubleCBFastDict.o
# 
# RooDoubleCBFastDict.o : RooDoubleCBFastDict.cxx
# 	g++ -c RooDoubleCBFastDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)
# 
# RooDoubleCBFastDict.cxx : RooDoubleCBFast.h LinkDef.h
# 	rootcint -f RooDoubleCBFastDict.cxx -c RooDoubleCBFast.h LinkDef.h
# 
# RooDoubleCBFast.o : RooDoubleCBFast.cc RooDoubleCBFast.h
# 	g++ -c RooDoubleCBFast.cc -fPIC $(CFLAGS) $(LDFLAGS)  -I./vdt			
# 