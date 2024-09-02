file="_reactor.c"

gcc -O2 "$file" -o "out_reactor" -lm

echo "Running code"

./out_reactor
