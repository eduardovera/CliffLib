TEMPLATE = app
CONFIG += console c++14 debug
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS+= -fopenmp -lpthread
QMAKE_LFLAGS +=  -fopenmp

HEADERS += \
    metric.h \
    multivector.hpp


SOURCES += main.cpp
LIBS += -L/usr/local/lib -ldlib -lpthread -fopenmp
