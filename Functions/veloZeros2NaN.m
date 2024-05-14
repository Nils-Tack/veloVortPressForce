function veloData = veloZeros2NaN(veloData)
veloData.uVelo(veloData.uVelo==0) = NaN;
veloData.vVelo(veloData.vVelo==0) = NaN;
veloData.UVspeed(veloData.UVspeed==0) = NaN;
end