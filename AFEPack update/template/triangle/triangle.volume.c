double triangle_volume(double ** v)
{
	return .5*((v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]));
};