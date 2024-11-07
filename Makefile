NUM_OF_TEST = 4
INPUT_PATH = "inputs/"
test_vectors:
	num=1; while [[ $$num -le $(NUM_OF_TEST) ]]; do\
		echo $$num;\
		mkdir -p inputs/test_$$num;\
		echo $$(INPUT_PATH);\
	 	python3 ./scripts/input_gen_mini_project.py -r 4 -c 16 -x 16 -y 16 -t "$(INPUT_PATH)/test_$$num/test_$$num";\
		((num = num+1));\
	done

mystery_test_vectors:
	num=1; while [[ $$num -le $(NUM_OF_TEST) ]]; do\
		echo $$num;\
		mkdir -p inputs/mystery_test_$$num;\
		echo $$(INPUT_PATH);\
	 	python3 ./scripts/input_gen_mini_project.py -r 4 -c 16 -x 16 -y 16 -t "$(INPUT_PATH)/mystery_test_$$num/mystery_test_$$num";\
		((num = num+1));\
	done
