def safe_sum(array=[]):
    """tail recursive safe summing of an array of elements if thats even possible in python
    can't belive numpy is not doing this by default..."""
    num = len(array) # number of elements in the array
    if num == 0:
        return float('NaN')
    elif num == 1:
        return array[0]
    elif num % 2 == 1:
        return safe_sum([ array[i] if i == num-1 else array[i]+array[i+1] for i in range(0,num,2) ])
    else:
        return safe_sum([array[i] + array[i + 1] for i in range(0, num, 2)])
        
def ordered_last_search(array, predicate):
    """find the index of the last element which satisfies the predicate in an ordered array"""
    num = len(array)  # number of elements in the array
    num_half = int(num/2)
    if num == 0:
        return - 1
    elif num == 1:
        return 0
    elif predicate(array[num_half]):
        return num_half + ordered_last_search(array[num_half:], predicate)
    else:
        return ordered_last_search(array[:num_half], predicate)
