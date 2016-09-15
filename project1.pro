TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

unix {
    LIBS += -lblas -llapack -larmadillo
} else {
    INCLUDEPATH += C:/Armadillo/include/
    LIBS += -LC:/Armadillo/examples/lib_win64/ -lblas_win64_MT
    LIBS += -LC:/Armadillo/examples/lib_win64/ -llapack_win64_MT
    LIBS += -lblas_win64_MT -llapack_win64_MT -larmadillo
}
SOURCES += main.cpp
