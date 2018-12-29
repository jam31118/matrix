SRC_DIR = ./src/
TEST_DIR = ./test/
TEST_SUBDIRS := $(wildcard $(TEST_DIR)/*/.)

.PHONY : all test $(TEST_SUBDIRS)
all: objects

test: $(TEST_SUBDIRS)

objects:
	$(MAKE) -C $(SRC_DIR)

install:
	$(MAKE) install -C $(SRC_DIR)

$(TEST_SUBDIRS):
	$(MAKE) -C $@

clean:
	$(MAKE) clean -C $(SRC_DIR)

clean-test:
	for dir in $(TEST_SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done

