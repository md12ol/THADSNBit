cmake_minimum_required(VERSION 3.16)
project(THADSNBit)
#project(THADSNBit8)
#project(THADSNBit16)
#project(THADSNBit24)
#project(THADSNBit32)

set(CMAKE_CXX_STANDARD 17)

add_executable(THADSNBit
        Epi/main.cpp
        Epi/setu.h
        Epi/setu.cpp
        Epi/stat.h
        Epi/stat.cpp
        BitSprayer/bitspray.cpp
        BitSprayer/bitspray.h)
#add_executable(THADSNBit8
#        Epi/main.cpp
#        Epi/setu.h
#        Epi/setu.cpp
#        Epi/stat.h
#        Epi/stat.cpp
#        BitSprayer/bitspray.cpp
#        BitSprayer/bitspray.h)
#add_executable(THADSNBit16
#        Epi/main.cpp
#        Epi/setu.h
#        Epi/setu.cpp
#        Epi/stat.h
#        Epi/stat.cpp
#        BitSprayer/bitspray.cpp
#        BitSprayer/bitspray.h)
#add_executable(THADSNBit24
#        Epi/main.cpp
#        Epi/setu.h
#        Epi/setu.cpp
#        Epi/stat.h
#        Epi/stat.cpp
#        BitSprayer/bitspray.cpp
#        BitSprayer/bitspray.h)
#add_executable(THADSNBit32
#        Epi/main.cpp
#        Epi/setu.h
#        Epi/setu.cpp
#        Epi/stat.h
#        Epi/stat.cpp
#        BitSprayer/bitspray.cpp
#        BitSprayer/bitspray.h)

#add_executable(TEST
#        Epi/test.cpp
#        Epi/setu.h
#        Epi/setu.cpp)

set_property(TARGET THADSNBit PROPERTY WORKING_DIRECTORY ./)
#set_property(TARGET THADSNBit8 PROPERTY WORKING_DIRECTORY ./)
#set_property(TARGET THADSNBit16 PROPERTY WORKING_DIRECTORY ./)
#set_property(TARGET THADSNBit24 PROPERTY WORKING_DIRECTORY ./)
#set_property(TARGET THADSNBit32 PROPERTY WORKING_DIRECTORY ./)