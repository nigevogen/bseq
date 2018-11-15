def fasta_formatted_string(name, sequence, description=None, line_width=None):
    string = '>' + name
    if description:
        string += ' ' + description
    string += '\n'
    if line_width:
        last = 0
        for i in range(line_width, len(sequence), line_width):
            string += sequence[i-line_width:i] + '\n'
            last = i
        string += sequence[last:]
        return string + '\n'
    string += sequence
    return string + '\n'