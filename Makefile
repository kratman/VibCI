

### Compiler settings

CXX=g++
CXXFLAGS=-static -fopenmp -O3
DEVFLAGS=-g -Wall
LDFLAGS=-I./src/ -I/usr/include/eigen3/

### Compile rules for users and devs

install:	title lovcibin compdone

clean:	title delbin compdone

### Rules for building various parts of the code

lovcibin:	
	@echo ""; \
	echo "### Compiling the LOVCI binary ###"
	$(CXX) $(CXXFLAGS) ./src/LOVCI.cpp -o lovci $(LDFLAGS)

compdone:	
	@echo ""; \
	echo "Done."; \
	echo ""

title:	
	@echo ""; \
	echo "#########################################################"; \
	echo "#                                                       #"; \
	echo "# Ladder Operator Vibrational Configuration Interaction #"; \
	echo "#                                                       #"; \
	echo "#########################################################"; \
	echo ""

delbin:	
	@echo ""; \
	echo '     ___'; \
	echo '    |_  |'; \
	echo '      \ \'; \
	echo '      |\ \'; \
	echo '      | \ \'; \
	echo '      \  \ \'; \
	echo '       \  \ \'; \
	echo '        \  \ \       <wrrr vroooom wrrr> '; \
	echo '         \__\ \________'; \
	echo '             |_________\'; \
	echo '             |__________|  ..,  ,.,. .,.,, ,..'; \
	echo ""; \
	echo ""; \
	echo "Removing executables..."; \
	rm -f lovci
