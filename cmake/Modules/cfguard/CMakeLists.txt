cmake_minimum_required(VERSION 3.15)

project(cfguard VERSION 0.1.0 LANGUAGES C)

add_library(cfguard STATIC mingw_cfguard_support.c)
