import random


def random_id(length: int = 10) -> str:
    chars = "23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"
    return "".join(random.choice(chars) for _ in range(length))
