from scipy import constants

prefactor = 1/( 4.184/(constants.value('elementary charge')**2 * 1e7) *
               constants.epsilon_0 / constants.Avogadro) # 1/epsilon_0