.PHONY: all test clean without_optimizations

all:
	stack build --test --no-run-tests --copy-bins --local-bin-path .

without_optimizations:
	stack build --fast --test --no-run-tests --copy-bins --local-bin-path .

test:
	stack build --fast --test

clean:
	stack clean
