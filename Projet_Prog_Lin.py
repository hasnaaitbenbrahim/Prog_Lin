import heapq

# This function, identity, creates an identity matrix with a specified number of rows and columns. 

def identity(numRows, numCols, val=1, rowStart=0):
    return [[(val if i == j else 0) for j in range(numCols)]
               for i in range(rowStart, numRows)]



#This function, standardForm, transforms a linear program from its general form to its standard form. 

def standardForm(cost, greaterThans=[], gtThreshold=[], lessThans=[], ltThreshold=[],
                equalities=[], eqThreshold=[], maximization=True):
    newVars = 0
    numRows = 0
    if gtThreshold != []:
        newVars += len(gtThreshold)
        numRows += len(gtThreshold)
    if ltThreshold != []:
        newVars += len(ltThreshold)
        numRows += len(ltThreshold)
    if eqThreshold != []:
        numRows += len(eqThreshold)

    if not maximization:
        cost = [-x for x in cost]
    if newVars == 0:
        return cost, equalities, eqThreshold

    newCost = list(cost) + [0] * newVars

    constraints = []
    threshold = []
    oldConstraints = [(greaterThans, gtThreshold, -1), (lessThans, ltThreshold, 1),
                     (equalities, eqThreshold, 0)]
    offset = 0
    for constraintList, oldThreshold, coefficient in oldConstraints:
        constraints += [c + r for c, r in zip(constraintList,
         identity(numRows, newVars, coefficient, offset))]

        threshold += oldThreshold
        offset += len(oldThreshold)

    return newCost, constraints, threshold

#Cette fonction calcule le produit scalaire de deux vecteurs a et b. 
#Elle multiplie chaque élément correspondant de a et b, puis additionne tous les résultats pour obtenir une seule valeur.


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


#Cette fonction extrait la colonne j de la matrice A. Pour cela, 
#elle parcourt chaque ligne de A et prend l'élément situé à l'indice j de chaque ligne, 
#puis les rassemble dans une liste pour former la colonne.


def column(A, j):
    return [row[j] for row in A]

#Cette fonction transpose la matrice A, c'est-à-dire qu'elle échange
# les lignes et les colonnes de la matrice. Pour cela, elle utilise
# la fonction column pour extraire chaque colonne de A, puis les rassemble pour former les lignes de la matrice transposée.


def transpose(A):
    return [column(A, j) for j in range(len(A[0]))]

#Cette fonction vérifie si une colonne donnée est une colonne pivot dans le tableau du simplexe.
# Une colonne pivot est une colonne qui contient un seul élément égal à 1 et tous les autres éléments sont nuls.
# De plus, la somme de tous les éléments de la colonne doit être égale à 1,
# indiquant ainsi qu'elle correspond à une variable de base dans le tableau.


def isPivotCol(col):
    return (len([c for c in col if c == 0]) == len(col) - 1) and sum(col) == 1


#Cette fonction trouve la valeur de la variable associée à une colonne pivot dans le tableau du simplexe. 
#Elle recherche l'indice de ligne où se trouve le seul élément égal à 1 dans la colonne, 
#puis retourne la valeur située à la dernière colonne de cette ligne, ce qui représente 
#la valeur de la variable dans la solution du problème d'optimisation linéaire.

def variableValueForPivotColumn(tableau, column):
    pivotRow = [i for (i, x) in enumerate(column) if x == 1][0]
    return tableau[pivotRow][-1]


#Cette fonction, initialTableau, crée le tableau initial pour le problème 
#de programmation linéaire dans le cadre de la méthode du simplexe.

def initialTableau(c, A, b):
    tableau = [row[:] + [x] for row, x in zip(A, b)]
    tableau.append([ci for ci in c] + [0])
    return tableau


#Cette fonction, primalSolution, calcule la solution primaire à partir du tableau du simplexe
def primalSolution(tableau):
   
    columns = transpose(tableau)
    indices = [j for j, col in enumerate(columns[:-1]) if isPivotCol(col)]
    return [(colIndex, variableValueForPivotColumn(tableau, columns[colIndex]))
            for colIndex in indices]


#Cette fonction, objectiveValue, calcule la valeur de la fonction objectif à partir du tableau du simplexe. 
def objectiveValue(tableau):
    return -(tableau[-1][-1])


#Cette fonction, canImprove, vérifie si une amélioration de la solution est possible à partir du tableau du simplexe actuel.
def canImprove(tableau):
    lastRow = tableau[-1]
    return any(x > 0 for x in lastRow[:-1])


#Cette fonction, moreThanOneMin, vérifie s'il existe plus d'un minimum dans 
#ne liste de tuples selon le deuxième élément de chaque tuple. 

def moreThanOneMin(L):
    if len(L) <= 1:
        return False

    x,y = heapq.nsmallest(2, L, key=lambda x: x[1])
    return x == y


def findPivotIndex(tableau):
    # cherche l'index minimum positif de la dernière ligne
    column_choices = [(i,x) for (i,x) in enumerate(tableau[-1][:-1]) if x > 0]
    column = min(column_choices, key=lambda a: a[1])[0]

    # vérifier si non borné
    if all(row[column] <= 0 for row in tableau):
        raise Exception("Le problème de programmation linéaire n'a pas de limite inférieure ou supérieure définie.")


    # vérifier la dégénérescence : plus d'un minimiseur du quotient
    quotients = [(i, r[-1] / r[column])
        for i,r in enumerate(tableau[:-1]) if r[column] > 0]

    if moreThanOneMin(quotients):
        raise Exception('Le problème de programmation linéaire est dans un état dégénéré.')

    # choisir l'index de ligne minimisant le quotient
    row = min(quotients, key=lambda x: x[1])[0]

    return row, column


#Cette fonction, pivotAbout, effectue l'opération de pivotement 
#autour d'un pivot spécifique dans le tableau du simplexe. 

def pivotAbout(tableau, pivot):
    i,j = pivot

    pivotDenom = tableau[i][j]
    tableau[i] = [x / pivotDenom for x in tableau[i]]

    for k,row in enumerate(tableau):
        if k != i:
            pivotRowMultiple = [y * tableau[k][j] for y in tableau[i]]
            tableau[k] = [x - y for x,y in zip(tableau[k], pivotRowMultiple)]


#Cette fonction, simplex, implémente l'algorithme du simplexe
# pour résoudre un problème de programmation linéaire.
 
def simplex(c, A, b):
    tableau = initialTableau(c, A, b)
    print("Tableau initial:")
    for row in tableau:
        print(row)
    print()

    while canImprove(tableau):
        pivot = findPivotIndex(tableau)
        print("Prochain index de pivotage est=%d,%d \n" % pivot)
        pivotAbout(tableau, pivot)
        print("Tableau après pivotage:")
        for row in tableau:
            print(row)
        print()

    return tableau, primalSolution(tableau), objectiveValue(tableau)



if __name__ == "__main__":
    num_variables = int(input("insérer  le nbr de variables: "))
    c = []
    for i in range(num_variables):
        c.append(float(input("insérer le coefficient correspondant à la variable %d dans la fonction objectif: " % (i+1))))
    
    num_constraints = int(input("insérer le nombre de contraintes: "))
    A = []
    for i in range(num_constraints):
        constraint = []
        for j in range(num_variables):
            constraint.append(float(input("insérer2 le coefficient correspondant à la variable %d dans la fonction objectif " % (j+1, i+1))))
        A.append(constraint)
    
    b = []
    for i in range(num_constraints):
        b.append(float(input("insérer la valeur de la contrainte %d: " % (i+1))))

    t, s, v = simplex(c, A, b)
    print("Résultat  optimale:")
    print(s)
    print("Valeur de la fonction objectif:")
    print(v)

# RÉALISE PAR :
# -	OUMAIMA AIT BIHI 
# -	HASNA AIT BEN BRAHIM
# ENCADRÉ PAR : 
# -	M. a. Reha
#  
