FLAGS = -I../../include -lmatrix -L../../lib 
all: main

main: main.cc libmatrix.a
	$(CXX) -o $@ $< $(FLAGS)

libmatrix.a:
	$(MAKE) install -C ../../src

clean:
	$(MAKE) clean -C ../../src
	$(RM) -f main

