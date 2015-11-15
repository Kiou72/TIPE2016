from random import randrange
import random
import math


# Génération de la clé

def is_prime(n):
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False

    return True


def random_prime(i, j):
    prime = False
    while not prime:
        n = random.randrange(i, j)
        x = 2
        prime = True
        while prime and x <= int(math.sqrt(n)) + 1:
            prime = not ((n % x) == 0)
            x += 1
    return n


def randomprime(n):
    r = randrange(2, n)
    while not is_prime(r):
        r = randrange(2, n)
    return r


def generate_primes_nb(l):
    x = random_prime(math.pow(10, l - 1), math.pow(10, l) - 1)
    y = random_prime(math.pow(10, l - 1), math.pow(10, l) - 1)
    while x == y:
        y = random_prime(math.pow(10, l - 1), math.pow(10, l) - 1)
    return x, y


def puiss_mod(x, p, n):
    """x^p % n"""
    puis = x % n
    for k in range(p, 1, -1):
        puis = (puis * x) % n
    return puis


def pgcd(a, b):
    r = a % b
    if r == 0:
        return b
    else:
        return pgcd(b, r)


def ppcm(a, b):
    return a*b//pgcd(a, b)


def bezout(a, b):
    """ Calcule (u, v, p) tels que a*u + b*v = p et p = pgcd(a, b) """
    r1 = a
    r2 = b
    u1 = 1
    u2 = 0
    v1 = 0
    v2 = 1
    while r2 != 0:
        q = r1//r2
        r1, r2 = r2, r1 - q * r2
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2
    return u1, v1, r1


def inv_modulo(x, m):
    """ Calcule y dans [[0, m-1]] tel que x*y % abs(m) = 1 """
    (u, _, p) = bezout(x, m)
    return u % m


def generate_key():
    (p, q) = generate_primes_nb(2)
    n = p * q
    n2 = math.pow(n, 2)
    lpp = ppcm(p - 1, q - 1)
    g = randrange(1, n2)
    glamb = puiss_mod(g, lpp, n2)
    l = (glamb - 1) // n
    mu = inv_modulo(l, n)
    return (n, g), (lpp, int(mu))


(pub, priv) = generate_key()
print(pub)
print(priv)


# Encryptage du message

def rand_znet(n):
    prime = False
    k = 0
    while not prime:
        k = random.randrange(-n + 1, n - 1)
        if pgcd(k, n) == 1:
            prime = True
    return k


def encoder(m, pub):
    (n, g) = pub
    r = rand_znet(n)
    n2 = n ** 2
    gm = puiss_mod(g, m, n2)
    rn = puiss_mod(r, n, n2)
    return (gm * rn) % (n2)


c = encoder(73, pub)
print(c)


# Décryptage du message

def decoder(m, priv, pub):
    (n, g) = pub
    (lamb, mu) = priv
    cl = puiss_mod(c, lamb, n ** 2)
    l = (cl - 1) // n
    return (l * mu) % n


print(decoder(c, priv, pub))

#for k in range(10000):
#    (pub, priv) = generate_key()
 #   c = encoder(73, pub)
  #  if decoder(c, priv, pub) != 73:
   #     print("faux")