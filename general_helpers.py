import pandas
from collections import Counter


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def plot_counter_list(l):
    letter_counts = Counter(l)
    df = pandas.DataFrame.from_dict(letter_counts, orient='index')
    df.plot(kind='bar')


def search_in_list(count_list, key):
    d = {v[0]: v[1] for v in count_list}
    return d.get(key, 0)


def is_overlap(l1 , l2):
    return len([f for f in l1 if f in l2])!=0




