
import math
import numpy as np

class AlbersEqualArea:

    def __init__(self):

        # geographic constants
        self.a = 6378137.0
        f = 1.0/298.2572221010042
        self.e2 = 2*f - f**2
        self.e = math.sqrt(self.e2)
        phi_0 = math.radians(23.0)
        phi_1 = math.radians(29.5)
        phi_2 = math.radians(45.5)
        self.lambda_0 = math.radians(-96.0)

        # preliminaries
        m_1 = self._m(phi_1)
        m_2 = self._m(phi_2)
        q_0 = self._q(phi_0)
        q_1 = self._q(phi_1)
        q_2 = self._q(phi_2)

        # map constants
        self.n = (m_1**2 - m_2**2)/(q_2 - q_1)
        self.C = m_1**2 + self.n*q_1
        self.rho_0 = self.a*(self.C - self.n*q_0)**0.5/self.n

    def _m(self, phi):

        return math.cos(phi)/math.sqrt(1 - self.e2*math.pow(math.sin(phi),2))

    def _q(self, phi):

        return (1 - self.e2)*(math.sin(phi)/(1 - self.e2*math.pow(
            math.sin(phi),2)) - (1.0/(2*self.e))*math.log((1-self.e*math.sin(
            phi))/(1 + self.e*math.sin(phi))))

    def project(self, lat, lon):

        lat = math.radians(lat)
        lon = math.radians(lon)

        q = self._q(lat)

        rho = self.a*(self.C - self.n*q)**0.5/self.n
        theta = self.n*(lon - self.lambda_0)

        x = rho*math.sin(theta)
        y = self.rho_0 - rho*math.cos(theta)

        return x,y

if __name__ == '__main__':

    test = AlbersEqualArea()

    #lat, lon = 39.2213, -120.0019
    lat1, lon1 = 40.0, -121.0
    lat2, lon2 = 39.0, -120.0
    test = AlbersEqualArea()
    x1, y1 = test.project(lat1, lon1)
    x2, y2 = test.project(lat2, lon2)
    print((np.abs(x2 - x1)*np.abs(y2 - y1))/1000000)
