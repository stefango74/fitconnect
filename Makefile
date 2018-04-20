CODE_DIR = src

.PHONY: code

code:
	$(MAKE) -C $(CODE_DIR)

clean:
	$(MAKE) -C $(CODE_DIR) clean
