CC=g++
BIN=bin

all: uavsar_calib uavsar_geocode

uavsar_calib: uavsar_calib.cpp load_ann.h math_uavsar.h optionparser.h
		$(CC) uavsar_calib.cpp -o uavsar_calib

uavsar_geocode: uavsar_geocode.cpp
		$(CC) uavsar_geocode.cpp -o uavsar_geocode

install:
	if [ ! -d $(BIN) ]; then mkdir $(BIN); fi; mv uavsar_calib $(BIN); mv uavsar_geocode $(BIN);

clean:
	rm -rf *.h.gch *~ uavsar_calib uavsar_geocode
