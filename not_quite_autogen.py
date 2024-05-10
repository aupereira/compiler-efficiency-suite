import json
from random import random, shuffle
import subprocess
import sys
import os


def get_bash_var(CMD):
    process = subprocess.Popen(CMD, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()

    if error:
        print(f"Error: {error}")
        sys.exit(1)
    else:
        return output.decode()


def load_json():
    with open("config.json") as f:
        data = json.load(f)
        p = data['profiler_cmd']
        r = data['runs_per_config']
        l = data['languages']
        b = data['benchmarks']
        return p, r, l, b


def gen_mk_dirs(BASE_DIR):
    output = ""

    if sys.platform == "linux":
        output += f"rm -rf {BASE_DIR}\n"

    for i in BENCHMARKS:
        for j in LANGUAGES:
            for k in j['compilers']:
                dir = os.path.join(BASE_DIR, f"{i['name']}_{k['name']}")
                output += f"mkdir --parents {dir}\n"

    return output


def gen_path(BASE_ENV, compiler, rand_len):
    source = ""
    size = ENV_SIZE
    
    if 'extra_env' in compiler:
        source += f"export {compiler['extra_env'].replace('$USER', USER)}\n"
        size += len(compiler['extra_env'])

    if 'path' in compiler:
        source += f"export PATH={compiler['path'].replace('$USER', USER)}:{BASE_ENV}"
        size += len(compiler['path'])

    if rand_len:
        source += "\"" + (" " * int(random() * (3584 - (size)))) + "\""

    source += "\n"

    return source


def gen_benchmark_compile():
    output = ""

    for i in BENCHMARKS:
        for j in LANGUAGES:
            for k in j["compilers"]:
                source_dir = os.path.join(SRC_DIR, i['name'], f"{i['name']}_{j['extension']}")
                build_dir = os.path.join(BUILD_DIR, f"{i['name']}_{k['name']}")
                
                output += f"echo Building {i['name']}_{k['name']}\n"

                output += gen_path(BASE_ENV, k, False)
                
                if 'build_cmd' in k:
                    output += f"cd {source_dir}\n"

                    compilers_with_args = [ind['name'] for ind in i['compiler_args']]

                    output += k['build_cmd'].replace("{output}", f"{i['name']}_{j['extension']}")
                    if k['name'] in compilers_with_args:
                        output += f" {i['compiler_args'][compilers_with_args.index(k['name'])]['args']}"

                    output += "\n"

                    if j['language'] in ['Java']:
                        exec_path = os.path.join(source_dir, "Main.class")
                    elif j['language'] in ['C#']:
                        exec_path = os.path.join(source_dir, k['output_dir']) + "/*"
                    elif 'output_dir' in k:
                        exec_path = os.path.join(source_dir, k['output_dir'], f"{i['name']}_{j['extension']}")
                    else:
                        exec_path = os.path.join(source_dir, f"{i['name']}_{j['extension']}")

                    output += f"mv {exec_path} {build_dir}\n"

                    if j['language'] not in ['Java']:
                        output += f"cd {build_dir}\n"
                        old_name = f"{i['name']}_{j['extension']}"
                        new_name = f"{i['name']}_{k['name']}"
                        output += f"mv {old_name} {new_name}\n"
                
                else:
                    script_path = os.path.join(source_dir, f"main.{j['extension']}")
                    output += f"cp {script_path} {build_dir}\n"
                    
                    output += f"cd {build_dir}\n"
                    old_name = f"main.{j['extension']}"
                    new_name = f"{i['name']}_{k['name']}.{j['extension']}"
                    output += f"mv {old_name} {new_name}\n"

                if 'extra_env' in k:
                    output += f"unset {k['extra_env'][:k['extra_env'].index('=')]}\n"

                output += "\n"

    return output


def gen_benchmark_run():
    outputs = []

    for iter in range(int(RUNS_PER_CONFIG)):
        for i in BENCHMARKS:
            output_buffer = []
            for j in LANGUAGES:
                for k in j['compilers']:
                    for config in range(len(i['test_args'])):
                        output = f"echo $(date)\n"
                        
                        benchmark = f"{i['name']}_{k['name']}"

                        output += gen_path(BASE_ENV, k, True)

                        output += "cd " + os.path.join(BUILD_DIR, f"{i['name']}_{k['name']}") + "\n"

                        # Hacky workaround for overriding the configs of slower languages.
                        if k['name'] in ['pypy', 'cpython', 'ruby', 'jruby'] and i['name'] == 'eigenvalue':
                            test_args = ['32 3 144 1', '32 3 64 12', '64 3 3 1', '64 3 2 12', '96 3 1 1', '96 3 1 12'][i['test_args'].index(i['test_args'][config])]
                        elif k['name'] in ['pypy', 'cpython', 'ruby', 'jruby', 'truffleruby'] and i['name'] == 'sha512':
                            test_args = i['test_args'][config].replace('test.txt', 'test_small.txt')
                        else: test_args = i['test_args'][config]

                        output += f"echo Running iteration {iter} of {i['name']}_{k['name']} with arguments {test_args}\n"

                        output += PROFILER_CMD.replace("{result_dir}", os.path.join(RESULTS_DIR, f"{i['name']}_{k['name']}", f"{config}_{iter}.txt"))

                        if 'run_cmd' in k:
                            output += f" {k['run_cmd'].replace('{file}', benchmark + '.' + j['extension'])} {test_args}"
                        else:
                            output += f" ./{benchmark} {test_args}"

                        output += "\n"
                        
                        if 'extra_env' in k:
                            output += f"unset {k['extra_env'][:k['extra_env'].index('=')]}\n"

                        output += "\n"
                        output_buffer.append(output)
            shuffle(output_buffer)
            outputs.append(''.join(output_buffer))
    
    return ''.join(outputs)


def main():
    if sys.platform == "linux":
        build_file = "build.sh"
        benchmark_file = "benchmark.sh"

    with open(build_file, "w") as f:
        f.write(gen_mk_dirs(BUILD_DIR))
        f.write("\n")
        f.write(gen_benchmark_compile())
        f.write(f"export PATH={BASE_ENV}")

    with open(benchmark_file, "w") as f:
        f.write(gen_mk_dirs(RESULTS_DIR))
        f.write("\n")
        f.write(gen_benchmark_run())
        f.write(f"export PATH={BASE_ENV}")

BASE_ENV = get_bash_var("echo $PATH").strip()
ENV_SIZE = len(get_bash_var("printenv").strip())
USER = get_bash_var("echo $USER").strip()

PROFILER_CMD, RUNS_PER_CONFIG, LANGUAGES, BENCHMARKS = load_json()

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
SRC_DIR = os.path.join(CURRENT_DIR, "src")
BUILD_DIR = os.path.join(CURRENT_DIR, "build")
RESULTS_DIR = os.path.join(CURRENT_DIR, "results")


if __name__ == "__main__":
    if sys.platform not in ["linux"]:
        print("This script only supports Linux")
        sys.exit(1)
    main()
