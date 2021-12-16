from Field import sum

def simple(fine, coarse):
    """Performs a simple contraction where every other value is deleted
    
    Parameters
    ----------
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    [x_fine, y_fine, dim] = fine.shape()

    # get slice indices
    xSlice = range(0,x_fine,2)
    ySlice = range(0,y_fine,2)

    # copy over
    iNew = 0
    for i in xSlice:
        jNew = 0
        for j in ySlice:
            if dim > 0:
                for k in range(dim):
                    coarse[iNew, jNew, k] = fine[i,j,k]
                coarse[iNew, jNew] = fine[i,j]
            jNew += 1
        iNew += 1


def sum4way(fine, coarse):
    """Contracts the Field by summing 4 values into 1
    
    Parameters
    ----------
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    # if not isinstance(fine, Field) or not isinstance(coarse, Field):
        # raise TypeError('Fine or coarse field is not a field')
    
    # Check that fine grid is divisible by 2 in both dims
    x_fine = fine.size()[0]
    y_fine = fine.size()[1]
    if (x_fine % 2) == 1 or (y_fine % 2) == 1:
        raise ValueError('Fine field dimensions do not allow for 4 way sum')

    # Check that dimensions of coarse grid are half of fine grid
    x_coarse = coarse.size()[0]
    y_coarse = coarse.size()[1]
    # if (x_fine / 2) != x_coarse or (y_fine / 2) != y_coarse:
    #     raise ValueError('Coarse grid size different from expected reduction from fine grid')
    
    shape = coarse.shape()
    dim = shape[2]

    ic = 0
    for i in range(0,x_fine,2):
        jc = 0
        for j in range(0,y_fine,2):
            for k in range(dim):
                coarse[ic, jc, k] = sum(fine[i:i+2,j:j+2,k])
            jc += 1
        ic += 1


def conservative4way(fine, coarse, weights=None):
    """Contracts the Field by averaging 4 values into 1 with weighting terms
    
    Parameters
    ----------
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    if weights is None:
        sum4way(fine, coarse)
        coarse.scale(0.25)
        return

    # Combines fine grid by summing over 4 blocks 
    # if not isinstance(fine, Field) or not isinstance(coarse, Field):
    #     raise TypeError('Fine or coarse field is not a field')
    
    # Check that fine grid is divisible by 2 in both dims
    x_fine = fine.size()[0]
    y_fine = fine.size()[1]
    if (x_fine % 2) == 1 or (y_fine % 2) == 1:
        raise ValueError('Fine field dimensions do not allow for 4 way sum')

    # Check that dimensions of coarse grid are half of fine grid
    x_coarse = coarse.size()[0]
    y_coarse = coarse.size()[1]
    # if (x_fine / 2) != x_coarse or (y_fine / 2) != y_coarse:
    #     raise ValueError('Coarse grid size different from expected reduction from fine grid')
    
    shape = coarse.shape()
    dim = shape[2]

    ic = 0
    for i in range(0,x_fine,2):
        jc = 0
        for j in range(0,y_fine,2):
            for k in range(dim):
                num = sum(fine[i:i+2,j:j+2,k] * weights[i:i+2,j:j+2])
                den = sum(weights[i:i+2,j:j+2])
                coarse[ic, jc, k] = num/den
            jc += 1
        ic += 1
