CFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags --glibs` -lRooFit -lRooFitCore -lMinuit

## https://root-forum.cern.ch/t/compiling-a-roofit-user-pdf-with-g/16266/3


libRooDoubleCBFast.so : libRooDoubleCBFast.so.1.0
	ln -sf libRooDoubleCBFast.so.1.0 libRooDoubleCBFast.so

libRooDoubleCBFast.so.1.0 : RooDoubleCBFast.o RooDoubleCBFastDict.o
	gcc -shared -Wl,-soname,libRooDoubleCBFast.so.1 -o libRooDoubleCBFast.so.1.0 RooDoubleCBFast.o RooDoubleCBFastDict.o

RooDoubleCBFastDict.o : RooDoubleCBFastDict.cxx
	g++ -c RooDoubleCBFastDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)

RooDoubleCBFastDict.cxx : RooDoubleCBFast.h LinkDef.h
	rootcint -f RooDoubleCBFastDict.cxx -c RooDoubleCBFast.h LinkDef.h

RooDoubleCBFast.o : RooDoubleCBFast.cc RooDoubleCBFast.h
	g++ -c RooDoubleCBFast.cc -fPIC $(CFLAGS) $(LDFLAGS)  -I./vdt	
	