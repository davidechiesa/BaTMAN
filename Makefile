LIBDIR  = lib/
INCDIR  = include/
BINDIR = bin/

ROOTCERN = `root-config --glibs --cflags`
CFLAGS = -shared -fPIC

CC = g++

all: $(LIBDIR)MyLib $(LIBDIR)NucReaction $(BINDIR)JAGS_input $(BINDIR)jags2root $(BINDIR)PosteriorAnalysis 
    
$(LIBDIR)MyLib: source/MyLib.cc source/date.cc $(INCDIR)MyLib.h $(INCDIR)date.h
	mkdir -p $(LIBDIR)  
	$(CC) $(CFLAGS) -o $(LIBDIR)libMyLib.so $+ -I $(INCDIR)
    
$(LIBDIR)NucReaction: source/NucReaction.cc $(INCDIR)NucReaction.h
	$(CC) $(CFLAGS) -o $(LIBDIR)libNucReaction.so $+ -I $(INCDIR) -L $(LIBDIR) -lMyLib $(ROOTCERN)
    
$(BINDIR)JAGS_input: source/JAGS_input.cpp
	mkdir -p $(BINDIR)
	$(CC) -o $(BINDIR)JAGS_input $+ -I $(INCDIR) -L $(LIBDIR) -lNucReaction -lMyLib $(ROOTCERN)
	
$(BINDIR)jags2root: source/jags2root.cpp source/jags2root.cc $(INCDIR)jags2root.h
	$(CC) -o $(BINDIR)jags2root $+ $(ROOTCERN) -I $(INCDIR) -L $(LIBDIR) -lMyLib 
	
$(BINDIR)PosteriorAnalysis: source/PosteriorAnalysis.cpp source/PosteriorAnalysis.cc $(INCDIR)/PosteriorAnalysis.h
	$(CC) -o $(BINDIR)PosteriorAnalysis $+ -I $(INCDIR) -L $(LIBDIR) -lNucReaction -lMyLib $(ROOTCERN)
	
clean: 
	rm $(LIBDIR)libMyLib.so
	rm $(LIBDIR)libNucReaction.so
	rm $(BINDIR)JAGS_input 
	rm $(BINDIR)jags2root
	rm $(BINDIR)PosteriorAnalysis



