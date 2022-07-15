for j in {0..10}
do
        taskset -c "$j" python3 main.py &
done