#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_FRI_PARAMETERS 16

void load_precomputed_inverses(const char *filename, uint64_t inverses[MAX_FRI_PARAMETERS]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[256];
    int domain;

    while (fgets(line, sizeof(line), file)) {
        // Look for the "FRI Domain" header
        if (sscanf(line, "FRI Domain %d Inverses:", &domain) == 1) {
            if (domain < 0 || domain >= MAX_FRI_PARAMETERS) {
                fprintf(stderr, "Invalid domain: %d\n", domain);
                exit(EXIT_FAILURE);
            }

            // Read only the **first inverse** after "FRI Domain" header
            if (fgets(line, sizeof(line), file)) {
                uint64_t value;
                if (sscanf(line, "%lx", &value) == 1) {
                    inverses[domain] = value;
                } else {
                    fprintf(stderr, "Failed to parse inverse for domain %d\n", domain);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    fclose(file);
}

int main() {
    uint64_t precomputed_inverses[MAX_FRI_PARAMETERS] = {0};

    const char *test_filename = "fri_inverses_256.txt";

    printf("Loading test file: %s\n", test_filename);
    load_precomputed_inverses(test_filename, precomputed_inverses);

    printf("\nLoaded precomputed inverses:\n");
    for (int i = 0; i < MAX_FRI_PARAMETERS; i++) {
        printf("Domain %d: %016lx\n", i, precomputed_inverses[i]);
    }

    return 0;
}