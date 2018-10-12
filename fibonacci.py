""" Can generate fibonacci sequence using user input as end limit. """


def end_integer_input():

    """

    Makes sure only positive integers above 1 are returned.
    Args: No args
    Returns: a positive whole integer user has entered only.

    """
    while True:
        response = input("Enter a number: ")
        try:
            end_value = int(response)
            if end_value >= 1:
                break
            else:
                print("Try again")
        except ValueError:
            print("Try again")
    return end_value


def generate_sequence(end_pos):

    """

    Generates a python list of subsequent fibonacci sequence until argument given.
    Args: a positive whole integer.
    Return: fibonacci sequence up until the integer position given as a list.

    """
    fib_sequence = [0]

    if end_pos >= 2:
        fib_sequence = [0,1]
        for number in range(end_pos-2):
            next_number_in_sequence = fib_sequence[number] + fib_sequence[number + 1]
            fib_sequence.append(next_number_in_sequence)
    print(fib_sequence)


def main():
    if __name__ == "__main__":
        end_pos = end_integer_input()
        generate_sequence(end_pos)

main()