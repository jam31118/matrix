LIB_DIR=../lib/


all: libmatrix.a

install: libmatrix.a
	mkdir -p $(LIB_DIR) && mv $^ $(LIB_DIR)

libmatrix.a: matrix.o
	$(AR) rcs $@ $^

matrix.o: matrix.cc
	$(CXX) -c $< -o $@

clean:
	$(RM) -f *.a *.o $(LIB_DIR)/*.a


