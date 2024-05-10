#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#define BUFFER_SIZE 4096

union Buffer
{
    char chars[BUFFER_SIZE];
    int ints[BUFFER_SIZE / sizeof(int)];
};

void genFile(char *path, long long size)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL)
    {
        printf("Error: Could not open file");
        exit(1);
    }

    union Buffer buffer;

    for (long long i = 0; i < size / BUFFER_SIZE; i++)
    {
        for (long long j = 0; j < BUFFER_SIZE / sizeof(int); j++)
        {
            buffer.ints[j] = rand();
        }
        fwrite(buffer.chars, 1, BUFFER_SIZE, f);
    }

    for (long long j = 0; j < (size % BUFFER_SIZE) / sizeof(int) + 1; j++)
    {
        buffer.ints[j] = rand();
    }

    fwrite(buffer.chars, 1, size % BUFFER_SIZE, f);

    fclose(f);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: ./filegen <File Path> <File Size [Bytes]> <Random Seed [True|False]>");
        exit(1);
    }

    if (strcmp(argv[3], "True") == 0 || strcmp(argv[3], "true") == 0)
    {
        srand(time(NULL));
    }

    genFile(argv[1], atoll(argv[2]));

    return 0;
}
