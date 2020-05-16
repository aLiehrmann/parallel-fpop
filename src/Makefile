SOURCES = Candidate.cpp Fpop.cpp Interval.cpp Ordered_list_of_intervals.cpp Quadratic.cpp rcpp_fpop.cpp main.cpp
OBJECTS = Candidate.o Fpop.o Interval.o Ordered_list_of_intervals.o Quadratic.o rcpp_fpop.o main.o
CXX = g++ -g
EXECUTABLE = fpop

FLAGS = -std=gnu++11 -DCPP -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I"/Library/Frameworks/R.framework/Versions/3.6/Resources/library/Rcpp/include" -I/usr/local/include -fPIC -Wall -O3 -fopenmp -Wall -Wextra -Wnon-virtual-dtor -Wpedantic -Wunused -Wlogical-op
FLAGS2 = -Wl,-headerpad_max_install_names -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation

all: $(EXECUTABLE)

main.o: main.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

rcpp_fpop.o: rcpp_fpop.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

Candidate.o: Candidate.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

Fpop.o: Fpop.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

Interval.o: Interval.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

Ordered_list_of_intervals.o: Ordered_list_of_intervals.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

Quadratic.o: Quadratic.cpp
	$(CXX) $(FLAGS) -c $^ -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(FLAGS) $(FLAGS2) $^ -o $@

clean:
	rm $(OBJECTS)
	rm fpop
