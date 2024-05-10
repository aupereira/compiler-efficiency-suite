import os
import csv
import re
import shutil
import statistics


ROOT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')
RESULTS_DIR = os.path.join(ROOT_DIR, 'results')
RESULTS_PARSED_DIR = os.path.join(ROOT_DIR, 'results_parsed')


def parse_directory(RESULTS_DIR):
    files = {}

    for dirpath, dirnames, filenames in os.walk(RESULTS_DIR):
        for filename in filenames:
            if filename.endswith('.txt'):
                path = os.path.join(dirpath, filename)
                parts = path.split(os.path.sep)
                benchmark = parts[-2].split('_')[0]
                config = parts[-1].split('_')[0]
                compiler = parts[-2].split('_')[1]
                run = parts[-1].split('_')[1].split('.')[0]

                with open(path, 'r') as file:
                    content = file.read().strip()
                    joules = re.search(r'([\d,]+\.\d+) Joules', content)
                    seconds = re.search(r'([\d,]+\.\d+) seconds', content)
                    if joules and seconds:
                        joules = float(joules.group(1).replace(',', ''))
                        seconds = float(seconds.group(1).replace(',', ''))

                        # Configuration Specific Adjustments
                        if benchmark == 'eigenvalue' and compiler in ['pypy', 'cpython', 'ruby', 'jruby']:
                            joules *= [(5120/144), (3072/64), (128/3), (96/2), (16/1), (12/1)][int(config)]
                            seconds *= [(5120/144), (3072/64), (128/3), (96/2), (16/1), (12/1)][int(config)]
                        elif benchmark == 'sha512' and compiler in ['pypy', 'cpython', 'ruby', 'jruby', 'truffleruby']:
                            joules *= 10_000_000_000 / 30_000_000
                            seconds *= 10_000_000_000 / 30_000_000

                if joules < 3 or seconds < 3:
                    print(f"Removing {benchmark} {config} {compiler} {run} due to low energy or runtime.")
                else:
                    file = f'{benchmark}_{config}'
                    if file not in files:
                        files[file] = {}
                    if compiler not in files[file]:
                        files[file][compiler] = [[], []]
                    files[file][compiler][0].append(joules)
                    files[file][compiler][1].append(seconds)
    
    return files


def run_statistics(files, total_rng, total_configs):
    for file in files:
    
        for compiler in files[file]:
            joules = statistics.mean(files[file][compiler][0])
            seconds = statistics.mean(files[file][compiler][1])
            watts = joules / seconds

            total_rng += (max(files[file][compiler][0])-min(files[file][compiler][0]))/min(files[file][compiler][0])
            total_configs += 1
            files[file][compiler] = [compiler, joules, seconds, watts]
        
        files[file] = [files[file][compiler] for compiler in files[file]]
    
    return files, total_rng, total_configs


def calculate_ratios(files):
    for file in files:
        min_joules = min(files[file], key=lambda x: x[1])[1]
        min_seconds = min(files[file], key=lambda x: x[2])[2]
        
        for i in range(len(files[file])):
            joules = files[file][i][1]
            seconds = files[file][i][2]
            
            energy_ratio = joules / min_joules
            time_ratio = seconds / min_seconds
            
            files[file][i].insert(1, energy_ratio)
            files[file][i].insert(2, time_ratio)
    
    return files


def sort(files):
    for file in files:
        files[file] = sorted(files[file], key=lambda x: x[3])
    return files


def format(files):
    for file in files:
        for i in range(len(files[file])):
            files[file][i][1] = f'{round(files[file][i][1], 2):.2f}'
            files[file][i][2] = f'{round(files[file][i][2], 2):.2f}'
            files[file][i][3] = f'{round(files[file][i][3], 2):.2f}'
            files[file][i][4] = f'{round(files[file][i][4], 2):.2f}'
            files[file][i][5] = f'{round(files[file][i][5], 2):.2f}'

    return files


def main():
    print(RESULTS_DIR)
    total_rng, total_runs = 0, 0
    
    files = parse_directory(RESULTS_DIR)
    files, total_rng, total_runs = run_statistics(files, total_rng, total_runs)
    files = calculate_ratios(files)
    files = sort(files)
    files = format(files)

    if os.path.exists(RESULTS_PARSED_DIR):
        shutil.rmtree(RESULTS_PARSED_DIR)

    os.mkdir(RESULTS_PARSED_DIR)

    for file in files:
        path = os.path.join(RESULTS_PARSED_DIR, f'{file}.csv')
        with open(path, 'w', newline='\n') as f:
                writer = csv.writer(f)
                writer.writerow(['compiler', 'energy_ratio', 'time_ratio', 'joules', 'seconds', 'watts'])
                for row in files[file]:
                    writer.writerow(row)


if __name__ == '__main__':
    main()
