attach("mf_secret_sharing.sage")

#this test FAILS!
#from_auditors is (hopefully partial) list of (auditor_id, coeff)
#we can actually restrict to public data with ``unique recon''
#this becomes a criterion that can be investigated further
def attack(from_auditors, pub):
    basis = pub[1]
    M = basis[0].parent()
    ring = M.base_ring()
    fsigma = [ring.zero() for ix in range(0,M.sturm_bound()+1)]
    for aud_pair in from_auditors:
        fsigma[aud_pair[0]] = aud_pair[1]
    cols = [[ring.zero() for ix in range(0,M.sturm_bound()+1)] for b in basis]
    for ix in range(0,len(basis)):
        for aud_pair in from_auditors:
            cols[ix][aud_pair[0]] = basis[ix].coefficients([aud_pair[0]])[0]
    cols.append(fsigma)
    coeff_matrix = Matrix(cols).transpose()
    #if we can reconstruct the mf with less information, we don't satisfy
    #the fundamental property of secret sharing as defined by Shamir
    assert coeff_matrix.rank() != M.dimension()

def test_recon_random(level, weight, prime, num_senders, r):
    ring = GF(prime)
    M = ModularForms(level, weight, base_ring = ring)
    pub = create_public_data(M, num_senders, r)
    votes = [ring(randint(0,1)) for ix in range(0,num_senders)]
    sum_votes = sum(votes)

    #senders send (partial) votes to auditors
    auditors = [ring.zero() for ix in range(0,M.sturm_bound()+1)]
    for sender in range(0,num_senders):
        for auditor in range(0,len(auditors)):
            #auditors collect (partial) votes
            auditors[auditor] += send_to_auditor(votes[sender], sender, auditor, pub)

    #auditors send (partial) votes to senders
    sender_recv=[[] for ix in range(0,num_senders)]
    for auditor in range(0,len(auditors)):
        for sender in range(0,num_senders):
            #need id because order is important, and comms might be async
            sender_recv[sender].append((auditor,auditors[auditor]))

    #senders reconstruct vote sum
    sum_votes_recon = [ring.zero() for ix in range(0,num_senders)]
    for sender in range(0,num_senders):
        sum_votes_recon[sender] = reconstruct(sender, sender_recv[sender], pub)
        assert sum_votes_recon[sender] == sum_votes, "Reconstructed vote sum is incorrect."

def runtests_random_r():
    for (lvl,wgt,pr) in [(23,2,7),(23,2,11),(29,2,7),(29,2,11),(127,2,7),(127,2,11),(127,2,13)]:
        M = ModularForms(lvl,wgt,base_ring=GF(pr))
        for num_senders in range(1,min(pr-1,M.dimension())+1):
            r = randint(100,200)
            try:
                test_recon_random(lvl,wgt,pr,num_senders,r)
            except ValueError as err:
                print("ValueError; {0}".format(err))
            except RuntimeError as err:
                #this is BAD!
                print("RuntimeError; {0}".format(err))
