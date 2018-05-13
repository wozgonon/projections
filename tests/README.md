# Integration tests

This is the [standard rust folder](https://doc.rust-lang.org/book/second-edition/ch11-03-test-organization.html) for integration tests.

## Download and render

I'm not quite sure were to put the following, whether this folder should only contain integration tests written in rust.
Perhaps they ought to be in the examples folder rather than here in the integration tests folder.

You need 'make' and 'curl' installed to run these tests, this generates a number of .html files in the target folder
which can be inspected in a web browser:

```
$ make -f run-download-and-render-tests.makefile clean all
```

