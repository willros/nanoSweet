#define NOB_IMPLEMENTATION
#include "nob.h"

// add this compiler option for shadow stack: -mshstk in gcc
// add -static to compiler option

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF(argc, argv);

    const char *program = nob_shift_args(&argc, &argv);

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "gcc");
    nob_cmd_append(&cmd, "-o", "nanodup");
    nob_cmd_append(&cmd, "nanodup.c", "thpool.c");
    nob_cmd_append(&cmd, "-lz", "-lpthread", "-lm", "-O3");
    if (!nob_cmd_run_sync(cmd)) return 1;

    cmd.count = 0;
    nob_cmd_append(&cmd, "gcc");
    nob_cmd_append(&cmd, "-o", "nanotrim");
    nob_cmd_append(&cmd, "nanotrim.c", "thpool.c");
    nob_cmd_append(&cmd, "-lz", "-lpthread", "-lm", "-O3");
    if (!nob_cmd_run_sync(cmd)) return 1;

    // nanomux
    cmd.count = 0;
    nob_cmd_append(&cmd, "c++", "-c", "edlib/src/edlib.cpp", "-o", "edlib.o", "-I", "edlib/include");
    if (!nob_cmd_run_sync(cmd)) return 1;
    
    cmd.count = 0;
    nob_cmd_append(&cmd, "cc", "-c", "nanomux.c", "-o", "nanomux.o", "-Iedlib/include", "-O3");
    if (!nob_cmd_run_sync(cmd)) return 1;
    
    cmd.count = 0;
    nob_cmd_append(&cmd, "cc", "-c", "thpool.c", "-o", "thpool.o", "-O3");
    if (!nob_cmd_run_sync(cmd)) return 1;
    
    cmd.count = 0;
    nob_cmd_append(&cmd, "c++", "nanomux.o", "thpool.o", "edlib.o", "-lz", "-lpthread", "-lm", "-o", "nanomux");
    if (!nob_cmd_run_sync(cmd)) return 1;
    
    return 0;
}
