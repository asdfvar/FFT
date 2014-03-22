CC=g++
OPT=-O3

CFILES = testfft.cpp fft.cpp
HEADERS = Complex.h Complex.cpp
OBS=testfft.o fft.o
EXE = testFFT

all: $(EXE)

$(EXE): $(OBS)
	$(CC) $(OPT) $(OBS) $(HEADERS) -o $(EXE)

testfft.o: testfft.cpp
	$(CC) $(OPT) testfft.cpp -c

fft.o: fft.cpp
	$(CC) $(OPT) fft.cpp -c

clean:
	rm $(EXE) $(OBS)
