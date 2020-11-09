from hyperopt import hp
from hyperopt import fmin, tpe, space_eval, Trials
from sklearn.metrics import accuracy_score

# define an objective function
def objective(args):
    case, val = args
    if case == 'case 1':
        return val
    else:
        return val ** 2


def some_weird_func(args):
    """
    3*13+20+30=89     opt params x=13, y=20, z=14
    :param x: [0-100] opt 20
    :param y: [0 - 100] opt 20
    :param z: [0 - 100] opt 1-14  if z>x
    :return:
    """
    x, y, z = args
    ret = 0
    if x>20:
        ret += -x
    else:
        ret += 3 * x

    if y>20:
        ret -= y
    elif y<5:
        ret -= y
    else:
        ret += y

    if z>x and z<15:
        ret += 30

    return -ret


def some_weird_func2(args):
    """
    3*13+20+30=89     opt params x=13, y=20, z=14
    :param x: [0-100] opt 20
    :param y: [0 - 100] opt 20
    :param z: [0 - 100] opt 1-14  if z>x
    :return:
    """
    x, y, z = args
    return x+y+z

def obj(args):

    return accuracy_score(y, predict)

if __name__ == '__main__':

    space = hp.choice('a',
                      [
                          ('case 1', 1 + hp.lognormal('c1', 0, 1)),
                          ('case 2', hp.uniform('c2', -10, 10))
                      ])
    list_space = [
        hp.uniform('a', 0, 50),
        hp.uniform('b', 0, 50),
     hp.uniform('c', 0, 50)]
    trials = Trials()

    best = fmin(some_weird_func, list_space, algo=tpe.suggest, max_evals=10, trials=trials)

    print(best)

    print(space_eval(list_space, best))
    print(trials.trials)


from xgboost import XGBClassifier
print(XGBClassifier())
"""
regularization = over-fit
max_depth [2, 30] >> too short under-fitting, too high over-fitting 
sub_sample = [0.1 - 1] >> regularization, which part of your data it'll take for model building. 1- all the data (not in use)
colsample_bytree = [0.1- 1] >>
colsample_bylevel = [0.1- 1] >> regularization, which part of columns (genomes) it'll take in each node. 1- all columns(not in use)
min_child_weight = [1 - 100] >> regularization,  limits the amount of min children in leaf, in regression XGBoost it's exactly the min number of cildren in leaf.
colsample_bynode
reg_alpha 
reg_lambda
n_estimators
max_depth = [10 - 1000] >> raising it makes it better till some point.
learning_rate = [0.01 - 0.8] >> important tuning - too small will never converge. too high not proper result.
gamma >> regularization
"""