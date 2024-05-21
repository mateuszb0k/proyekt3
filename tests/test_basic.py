from audioop import *

from build.Debug.cmake_example import *


freq = int(input("Podaj czestotliwosc: "))
print(sinus(freq))
print(cosinus(freq))
print(pila(freq))
print(prostokatny(freq))