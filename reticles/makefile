CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo
all: CNT_generator.exe CNT+Cu_generator.exe
CNT_generator.exe : CNT_generator_0.1.cpp
	$(CC) CNT_generator_0.1.cpp -o CNT_generator.exe $(AFLAGS) $(CFLAGS)

CNT+Cu_generator.exe : CNT+Cu_generator_0.1.cpp
	$(CC) CNT+Cu_generator_0.1.cpp -o CNT+Cu_generator.exe $(AFLAGS) $(CFLAGS)