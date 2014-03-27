CC = mpicc

BIN = numgens

all: $(BIN)

numgens: numgens.c
	$(CC) -o $@ -Wall $<

clean:
	rm -f $(BIN)
