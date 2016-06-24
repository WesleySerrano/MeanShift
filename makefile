COMPILER = g++
ANN_DIR = ~/Downloads/ann_1.1.2
ANN = -I $(ANN_DIR)/include -L $(ANN_DIR)/lib -lANN
	
main: main.cpp utils.h
	$(COMPILER) -g -o main main.cpp utils.h -lX11 -lpthread -lm $(ANN)