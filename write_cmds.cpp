#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>

unsigned number_of_digits(unsigned n)
{
    unsigned digits = 0;

    do
    {
        n /= 10;
        digits++;
    }
    while (n != 0);

    return digits;
}

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "no arg given\n";
        return 1;
    }

    unsigned n = std::atoi(argv[1]);
    unsigned n_digits = number_of_digits(n);

    char buf[10];
    char fmt[10];

    std::sprintf(fmt, "%%0%uu\n", n_digits);

    std::ofstream cmds("cmds");

    for(unsigned i = 1; i <= n; i++)
    {
        std::sprintf(buf, fmt, i);
        cmds << buf;
    }
}
