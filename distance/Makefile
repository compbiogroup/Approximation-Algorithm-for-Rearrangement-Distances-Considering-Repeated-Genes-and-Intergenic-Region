.PHONY: clean test debug all

all:
	cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
	cmake --build build
	mv build/dist .

debug:
	cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
	cmake --build build
	mv build/dist .

test:
	RC_PARAMS="max_size=100 verbose_progress=1" build/dist_test

clean:
	cd build && make clean
