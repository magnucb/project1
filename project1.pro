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

DISTFILES += \
    ../build-project1-Desktop_Qt_5_7_0_GCC_64bit-Debug/test.txt \
    ../build-project1-Desktop_Qt_5_7_0_GCC_64bit-Debug/test_copy.txt
