# Default target
all: run plot

# Step 1: Make execute.sh executable
prepare:
	chmod +x ./execute.sh

# Step 2: Run the solver + perf data generation
run: prepare
	./execute.sh

# Step 3: Plot using python
plot:
	python plot_metrics.py

# Clean: delete temp_dir
clean:
	rm -rf temp_dir

# Optional: also delete plots if needed
clean_all:
	rm -rf temp_dir plots
