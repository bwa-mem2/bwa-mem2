
ALL_TESTS=$(shell ./find_tests.py)
ALL_RESULTS=$(patsubst %, results/%_result, $(ALL_TESTS))

.PHONY: all
all: results $(ALL_RESULTS)
	@echo "============== TEST RESULTS ==============="
	python3 report.py

results:
	mkdir results

wiki_ts_200M_uint64:
	curl -L https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/SVN8PI | zstd -d > $@

osm_cellids_200M_uint64:
	curl -L https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/8FX9BV | zstd -d > $@

rmi: $(shell find ../src/) $(shell find ../rmi_lib/src/)
	cd .. && cargo build --release
	cp ../target/release/rmi .


results/%_result: results rmi wiki_ts_200M_uint64 osm_cellids_200M_uint64 FORCE
	@echo "Running test $*"
	@cd $* && make -s
	@cp $*/result results/$*_result
	@echo "Test $* finished."


.PHONY: clean FORCE
clean:
	rm -rf results rmi
	for test in $(ALL_TESTS) ; do make -C $$test clean ; done
FORCE: 
