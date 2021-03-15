genrule(
    name = "makeconverge",
    srcs = ["include", "main.cpp"],
    outs=["calculator"],
    cmd = "mpicxx -std=c++11 -o $@ -I$(SRCS)",
visibility=["//visibility:public"]
)
