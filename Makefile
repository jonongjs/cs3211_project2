CC = mpicc

BIN = project2

all: $(BIN)

numgens: numgens.c
	$(CC) -o $@ -Wall $<

project2: project2.c
	$(CC) -o $@ -Wall $<

clean:
	rm -f $(BIN)
