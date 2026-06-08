.PHONY: build release install check fmt test clean dist tag help

VERSION := $(shell cargo metadata --no-deps --format-version 1 2>/dev/null \
	| sed -n 's/.*"version":"\([^"]*\)".*/\1/p' | head -1)
TARGET  := $(shell rustc -vV | sed -n 's/host: //p')
ARCHIVE := preflightx-v$(VERSION)-$(TARGET)

build:
	cargo build

release:
	cargo build --release

install:
	cargo install --path . --locked

check:
	cargo clippy --all-targets -- -D warnings
	cargo fmt --check

fmt:
	cargo fmt

test:
	cargo test

# Local host-platform archive (binary + README + LICENSE); CI builds the full matrix on tag.
dist: release
	@mkdir -p dist/$(ARCHIVE)
	@cp target/release/preflightx dist/$(ARCHIVE)/
	@cp README.md LICENSE dist/$(ARCHIVE)/
	@cd dist && tar czf $(ARCHIVE).tar.gz $(ARCHIVE)
	@cd dist && shasum -a 256 $(ARCHIVE).tar.gz > $(ARCHIVE).tar.gz.sha256
	@rm -rf dist/$(ARCHIVE)
	@echo "Wrote dist/$(ARCHIVE).tar.gz (+ .sha256)"

tag:
	@echo "Tagging v$(VERSION) and pushing to trigger the release workflow..."
	git tag -a v$(VERSION) -m "Release v$(VERSION)"
	git push origin v$(VERSION)

clean:
	cargo clean
	rm -rf dist

help:
	@echo "Targets:"
	@echo "  build     Debug build"
	@echo "  release   Optimised release build"
	@echo "  install   cargo install to ~/.cargo/bin"
	@echo "  check     clippy + fmt --check"
	@echo "  fmt       Auto-format"
	@echo "  test      Run tests"
	@echo "  dist      Package host-platform release archive into dist/"
	@echo "  tag       Tag v$(VERSION) and push to trigger the CI release"
	@echo "  clean     cargo clean + remove dist/"
