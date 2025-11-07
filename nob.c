#define NOB_IMPLEMENTATION
#define NOB_STRIP_PREFIX
#include "nob.h"

Nob_Cmd cmd = {0};

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF(argc, argv);

    cmd_append(&cmd, "cc");
    cmd_append(&cmd, "-o", "nanodup");
    cmd_append(&cmd, "nanodup.c", "thpool.c");
    cmd_append(&cmd, "-lz", "-lpthread", "-O3");
    if (!cmd_run(&cmd)) return 1;

    cmd_append(&cmd, "cc");
    cmd_append(&cmd, "-o", "nanotrim");
    cmd_append(&cmd, "nanotrim.c", "thpool.c");
    cmd_append(&cmd, "-lz", "-lm", "-lpthread", "-O3");
    if (!cmd_run(&cmd)) return 1;


    cmd_append(&cmd, "cc");
    cmd_append(&cmd, "-o", "nanomux");
    cmd_append(&cmd, "nanomux.c", "thpool.c");
    cmd_append(&cmd, "-lz", "-lpthread", "-lm", "-O3");
    if (!cmd_run(&cmd)) return 1;
    return 0;
}
