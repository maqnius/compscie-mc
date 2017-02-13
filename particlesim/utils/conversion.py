from scipy import constants

prefactor = 1/( 4.184/(constants.value('elementary charge')**2 * 1e7) *
               constants.epsilon_0 / constants.Avogadro) # 1/epsilon_0

def beta_to_kelvin(beta):
    """

    Parameters
    ----------
    beta : float
        beta of Boltzman distribution in kcal/(mol*K)

    Returns
    -------
    T : float
        Temperature in Kelvin
    """
    return 1/(constants.R * 0.239006 * 1e-3 * beta)

def kelvin_to_beta(T):
    """
    Parameters
    ----------
    T : float
        Temperature in Kelvin

    Returns
    -------
    beta : float
        beta of Boltzman distribution in kcal/(mol*K)
        """
    return 1 / (constants.R * 0.239006 * 1e-3 * T)