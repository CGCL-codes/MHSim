set -e
make
rm main.o
#rm /usr/lib/libneuro.so
rm libneuro.so
g++ -shared -o libneuro.so NeuroSim/*.o neuro/*.o *.o
#cp libneuro.so /usr/lib/
echo "Completed!"
