import gzip


def open_text_maybe_gzip(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)
