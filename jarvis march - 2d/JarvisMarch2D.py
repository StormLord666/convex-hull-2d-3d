import sys
import numpy as np
import matplotlib.pyplot as plt

# Ελέγχουμε τις γωνίες των ΑΓΒ και ΑΒΓ και ποια είναι μεγαλύτερη. Αν η ΑΓΒ είναι μικρότερη, τότε κινούμαστε με τη φορά του ρολογιού (clockwise):
def cwTurn(p1, p2, p3):
	if (p3[1]-p1[1])*(p2[0]-p1[0]) <= (p2[1]-p1[1])*(p3[0]-p1[0]):
		return True
	return False

# Κύρια συνάρτηση τοο αλγορίθμου για τη λίστα σημείων S:
def JarvisMarch(S):
	plt.figure()  # Θέτουμε το figure για το plot μετά
	index = 0 
	n = len(S) # Θέτουμε πλήθος των σημείων στο S
	P = [None] * n # Θέτουμε λίστα σημείων που θα αποτελούν το convex hull
	l = np.where(S[:,0] == np.min(S[:,0])) # Θέτουμε l για έλεγχο αν είναι αληθές όταν το S είναι το μικρότερο δυνατό σημείο
	pointOnHull = S[l[0][0]] # Θέτουμε το "αριστερότερο "σημείο του S, που θα ανήκει στο convex hull
	i = 0
	#Σκοπός των παραπάνω είναι να θέσουμε το αριστερότερο σημείο. 
	while True:
		P[i] = pointOnHull  
		endpoint = S[0] # Θέτουμε αρχικό τελικό σημείο για ένα υποψήφιο άκρο στο convex hull
		for j in range(1,n):
			if (endpoint[0] == pointOnHull[0] and endpoint[1] == pointOnHull[1]) or not cwTurn(S[j],P[i],endpoint): #Ελέγχουμε για μια σπάνια περίπτωση και μπορεί να συμβεί μόνο όταν j == 1 και δεν έχει οριστεί ακόμα καλύτερο τελικό σημείο για τον βρόχο.
				endpoint = S[j] #Έπειτα από έλεγχο, τίθεται endpoint = S[j], δηλαδή βρίσκεται μεγαλύτερη στροφή αριστερά και ενημερώνεται το τελικό σημείο.
		i = i + 1
		pointOnHull = endpoint #Θέτουμε το σημείο αριστερότερα του αρχικού ως το νέο σημείο που μελετάμε.
		J = np.array([P[k] for k in range(n) if P[k] is not None])
		plt.clf()               # Εκκαθάριση του plot
		plt.plot(J[:,0],J[:,1], 'r-', picker=10)   # Σχεδιασμός γραμμών
		plt.plot(S[:,0],S[:,1],".k")              # Σχεδιασμός κουκκίδων
		plt.axis('off')         # Αφαιρούνται οι άξονες
		plt.show(block=False)   # Κλείσιμο του plot
		plt.pause(0.5)    # Μικρή παύση πριν το κλεισιμό του
		index += 1
		if endpoint[0] == P[0][0] and endpoint[1] == P[0][1]: #Ελέγχουμε αν επέστρεψε στο αρχικό σημείο ο hull. 
			break #Αν ναι, τότε η loop σταματά.
	for i in range(n):
		if P[-1] is None:
			del P[-1]
	P = np.array(P)
	
	#Ο ακόλουθος κώδικας αποσκοπεί στην ολοκλήρωση του plot. Δίχως αυτόν, το convex hull μένει μη ολοκληρωμένο, σταματώντας στο σημείο πριν κλείσει πλήρως.
	plt.clf()
	plt.plot(P[:,0],P[:,1], 'r-', picker=10)
	plt.plot([P[-1,0],P[0,0]],[P[-1,1],P[0,1]], 'r-', picker=10)
	plt.plot(S[:,0],S[:,1],".k")
	plt.axis('off')
	plt.show(block=False)
	plt.pause(0.5)
	return P

def main():
	try:
		N = int(sys.argv[1])
	except:
		N = int(input("Number of points: ")) #Ζητείται από τον χρήστη o ορισμός του Ν, το πλήθος δλδ των σημείων. 
  
	# Δημιουργούνται Ν τυχαία σημεία με συντεταγμένες στο διάστημα [0,200)x[0,200):
	P = np.array([(np.random.randint(0,200),np.random.randint(0,200)) for i in range(N)])
	L = JarvisMarch(P) #Εκτελείται η JarvisMarch
	
	# Καλούμε την matplotlib για να παραμείνει η εικόνα του τελειωμένου convex hull.
	plt.plot(L[:,0],L[:,1], 'b-', picker=5)
	plt.plot([L[-1,0],L[0,0]],[L[-1,1],L[0,1]], 'b-', picker=5)
	plt.plot(P[:,0],P[:,1],".r")
	plt.axis('off')
	plt.show()

if __name__ == '__main__':
	main()
